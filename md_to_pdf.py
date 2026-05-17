"""
md_to_pdf.py
=============

Markdown ファイルを PDF に変換するユーティリティ。

役割 : .md → HTML(日本語フォント対応 CSS 付き)→ ヘッドレス Edge で印刷 → PDF。
       pandoc 等の追加インストール不要(Python の markdown ライブラリのみ使用)。
使い方: python md_to_pdf.py <入力.md> [出力.pdf]
       出力を省略すると入力と同じ場所・同名の .pdf を生成する。
作成 : 2026-05-18
"""

import sys
import subprocess
import pathlib
import tempfile
import markdown

# Windows の Microsoft Edge(Chromium 版)
EDGE = r"C:\Program Files (x86)\Microsoft\Edge\Application\msedge.exe"

# 日本語対応・印刷向け CSS
CSS = """
@page { size: A4; margin: 16mm 14mm; }
body { font-family: "Yu Gothic","Meiryo","Noto Sans CJK JP",sans-serif;
       font-size: 10.5pt; line-height: 1.6; color: #111; }
h1 { font-size: 18pt; border-bottom: 2px solid #333; padding-bottom: 4px; }
h2 { font-size: 13.5pt; border-bottom: 1px solid #999; padding-bottom: 3px;
     margin-top: 1.3em; }
h3 { font-size: 11.5pt; color: #223a55; margin-top: 1.0em; }
table { border-collapse: collapse; width: 100%; margin: 0.6em 0; font-size: 8.7pt;
        table-layout: fixed; }
th, td { border: 1px solid #888; padding: 4px 6px; vertical-align: top;
         word-break: break-word; overflow-wrap: anywhere; }
th { background: #ececec; }
code { background: #f0f0f0; padding: 1px 3px; border-radius: 3px; font-size: 9pt; }
blockquote { border-left: 3px solid #bbb; margin: 0.6em 0; padding: 2px 10px;
             background: #f7f7f7; font-size: 9.5pt; }
hr { border: none; border-top: 1px solid #ccc; margin: 1.0em 0; }
"""


def main():
    if len(sys.argv) < 2:
        print("usage: python md_to_pdf.py <入力.md> [出力.pdf]")
        sys.exit(1)

    src = pathlib.Path(sys.argv[1]).resolve()
    if not src.exists():
        print(f"[エラー] 入力が見つかりません: {src}")
        sys.exit(1)
    pdf = (pathlib.Path(sys.argv[2]).resolve() if len(sys.argv) > 2
           else src.with_suffix(".pdf"))
    html_path = src.with_suffix(".html")

    # Markdown → HTML(表・コードブロック対応)
    text = src.read_text(encoding="utf-8")
    body = markdown.markdown(text, extensions=["tables", "fenced_code", "sane_lists"])
    html = ("<!doctype html><html lang='ja'><head><meta charset='utf-8'>"
            f"<style>{CSS}</style></head><body>{body}</body></html>")
    html_path.write_text(html, encoding="utf-8")

    # ヘッドレス Edge で HTML → PDF(既存 Edge と干渉しないよう一時プロファイル)
    with tempfile.TemporaryDirectory() as profile:
        subprocess.run(
            [EDGE, "--headless", "--disable-gpu", "--no-pdf-header-footer",
             f"--user-data-dir={profile}",
             f"--print-to-pdf={pdf}", html_path.as_uri()],
            check=True, timeout=120,
        )

    print(f"  HTML 中間ファイル: {html_path}")
    print(f"  -> PDF 生成完了  : {pdf}")


if __name__ == "__main__":
    main()
