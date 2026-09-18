"""Execute trusted tutorial Python examples and embed their actual outputs.

Run from the package directory with PYTHONPATH=src and DIFFPY_DATA_DIR set.
Shell setup commands are never executed by this renderer.
"""

from __future__ import annotations

import argparse
import contextlib
import html
import io
import re
import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from generate_tutorial_plots import embed_images


def render_outputs(html_path: Path) -> None:
    """Run examples in document order, sharing the tutorial's Python namespace."""
    source = html_path.read_text(encoding="utf-8")
    source = re.sub(
        r'\n<!-- output:example-\d+ -->.*?<!-- /output:example-\d+ -->',
        "", source, flags=re.DOTALL,
    )
    pattern = re.compile(
        r'<pre data-example="(example-\d+)"([^>]*)><code>(.*?)</code></pre>',
        re.DOTALL,
    )
    namespace = {"__name__": "__tutorial__"}

    def execute(match: re.Match[str]) -> str:
        example, attributes, code = match.groups()
        print(f"Executing {example}", flush=True)
        figure_match = re.search(r'data-figure="([^"]+)"', attributes)
        shown = 0

        def show(*args, **kwargs):
            nonlocal shown
            if figure_match is None:
                raise RuntimeError(f"{example} has no output figure path")
            numbers = plt.get_fignums()
            if len(numbers) != 1 or shown:
                raise RuntimeError(f"Expected one figure for {example}")
            figure_path = html_path.parent / figure_match[1]
            plt.figure(numbers[0]).savefig(
                figure_path, dpi=180, bbox_inches="tight", facecolor="white"
            )
            plt.close("all")
            shown += 1

        output = io.StringIO()
        original_show = plt.show
        plt.show = show
        try:
            with contextlib.redirect_stdout(output), warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always", RuntimeWarning)
                exec(compile(html.unescape(code), f"tutorial:{example}", "exec"), namespace)
        finally:
            plt.show = original_show
            plt.close("all")
        if figure_match and shown != 1:
            raise RuntimeError(f"No figure displayed by {example}")
        warning_text = "".join(
            f"{warning.category.__name__}: {warning.message}\n" for warning in caught
        )
        text = warning_text + output.getvalue()
        if not text.strip():
            text = "Figure output shown below.\n" if shown else "No text output.\n"
        rendered = (
            f'\n<!-- output:{example} -->\n'
            '<p class="output-label">Output</p>\n'
            f'<pre class="code-output"><code>{html.escape(text.rstrip())}</code></pre>\n'
            f'<!-- /output:{example} -->'
        )
        return match[0] + rendered

    updated, count = pattern.subn(execute, source)
    if count == 0:
        raise RuntimeError("No tutorial examples found")
    html_path.write_text(updated, encoding="utf-8")
    embed_images(html_path)
    print(f"Updated outputs for {count} Python examples.", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--html", type=Path, default=Path("docs/index.html"))
    render_outputs(parser.parse_args().html)
