#!/usr/bin/env python3
# ---- Point the moved 1.x site at its own folder ----
# Usage: python3 -I tools/relink_v1.py <docs/v1 folder> <old site root>
# Example: python3 -I tools/relink_v1.py docs/v1 https://jakebrownscombe.github.io/positionR
#
# The 1.x site was built for the root it names and is served from <root>/v1/.
# Its search index, sitemap, pkgdown.yml and some HTML attributes hold absolute
# URLs to the old root. Those URLs would open 2.0 pages. This rewrites each one
# to the same path under <root>/v1/. A URL already under /v1/ is left alone, so
# a second run changes nothing. The script prints the counts by kind, then
# exits 1 if any URL to the old root is left.

import pathlib
import re
import sys

KINDS = (".html", ".json", ".xml", ".yml")


# ---- Read the arguments ----
def read_args(argv):
    if len(argv) != 3:
        sys.exit("usage: python3 -I tools/relink_v1.py <docs/v1 folder> <old site root>")
    folder = pathlib.Path(argv[1])
    if not folder.is_dir():
        sys.exit(f"relink_v1: {folder} is not a folder")
    root = argv[2].rstrip("/")
    if not root.startswith("https://"):
        sys.exit(f"relink_v1: the site root must start with https://, got {root}")
    return folder, root


# ---- Find URLs to the old root ----
# The root must end a path segment, so positionR does not match positionR2.
def old_root_pattern(root):
    return re.compile(re.escape(root) + r"(?=[/\"'<>\s)]|$)(?!/v1(?:/|$))/?")


# ---- Rewrite the files ----
def rewrite(folder, root):
    pattern = old_root_pattern(root)
    new_root = root + "/v1/"
    counts = {kind: [0, 0] for kind in KINDS}
    for path in sorted(folder.rglob("*")):
        if not path.is_file() or path.suffix not in KINDS:
            continue
        text = path.read_text(encoding="utf-8")
        new_text, n = pattern.subn(new_root, text)
        if n:
            path.write_text(new_text, encoding="utf-8")
            counts[path.suffix][0] += 1
            counts[path.suffix][1] += n
    return counts


# ---- Count what is left ----
def count_left(folder, root):
    pattern = old_root_pattern(root)
    left = 0
    for path in folder.rglob("*"):
        if path.is_file() and path.suffix in KINDS:
            left += len(pattern.findall(path.read_text(encoding="utf-8")))
    return left


def main(argv):
    folder, root = read_args(argv)
    counts = rewrite(folder, root)
    for kind in KINDS:
        files, urls = counts[kind]
        print(f"relink_v1: {kind:5} {urls:4d} URLs rewritten in {files:3d} files")
    total_files = sum(c[0] for c in counts.values())
    total_urls = sum(c[1] for c in counts.values())
    print(f"relink_v1: total {total_urls} URLs rewritten in {total_files} files")
    left = count_left(folder, root)
    print(f"relink_v1: {left} URLs to {root}/ left outside /v1/")
    return 1 if left else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
