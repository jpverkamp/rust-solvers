#!/usr/bin/env python3
"""
strip_leading_lines.py

Recursively processes all text files in a directory and removes any leading
lines that are either empty or start with '#'.

Usage:
    python strip_leading_lines.py <directory> [--dry-run] [--ext .txt .md ...]

Arguments:
    directory       Path to the root directory to process
    --dry-run       Preview changes without modifying files
    --ext           File extensions to process (default: all files)

Examples:
    python strip_leading_lines.py ./my_docs
    python strip_leading_lines.py ./my_docs --dry-run
    python strip_leading_lines.py ./my_docs --ext .txt .md .py
"""

import argparse
import sys
from pathlib import Path


def strip_leading_lines(text: str) -> str:
    """Remove leading lines that are empty or start with '#'."""
    lines = text.splitlines(keepends=True)
    i = 0
    for line in lines:
        stripped = line.rstrip('\n').rstrip('\r')
        if stripped == '' or stripped.startswith('#'):
            i += 1
        else:
            break
    return ''.join(lines[i:])


def is_binary(path: Path) -> bool:
    """Quick check: read a chunk and look for null bytes."""
    try:
        chunk = path.read_bytes()[:8192]
        return b'\x00' in chunk
    except OSError:
        return True


def process_directory(root: Path, dry_run: bool, extensions: list[str]) -> None:
    files = sorted(root.rglob('*'))
    changed = skipped = unchanged = 0

    for path in files:
        if not path.is_file():
            continue

        # Filter by extension if specified
        if extensions and path.suffix.lower() not in extensions:
            continue

        # Skip binary files
        if is_binary(path):
            print(f"  [skip binary]  {path}")
            skipped += 1
            continue

        try:
            original = path.read_text(encoding='utf-8', errors='replace')
        except OSError as e:
            print(f"  [read error]   {path}: {e}")
            skipped += 1
            continue

        modified = strip_leading_lines(original)

        if modified == original:
            print(f"  [unchanged]    {path}")
            unchanged += 1
        else:
            removed = original.count('\n') - modified.count('\n')
            # Count leading lines removed
            original_lines = original.splitlines()
            modified_lines = modified.splitlines()
            lines_removed = len(original_lines) - len(modified_lines)
            label = "[dry-run]" if dry_run else "[modified]"
            print(f"  {label}  {path}  (-{lines_removed} leading line{'s' if lines_removed != 1 else ''})")
            if not dry_run:
                try:
                    path.write_text(modified, encoding='utf-8')
                except OSError as e:
                    print(f"  [write error]  {path}: {e}")
                    skipped += 1
                    continue
            changed += 1

    print()
    print(f"Done. {changed} file(s) {'would be ' if dry_run else ''}modified, "
          f"{unchanged} unchanged, {skipped} skipped.")


def main():
    parser = argparse.ArgumentParser(
        description="Strip leading empty/comment lines from text files recursively."
    )
    parser.add_argument('directory', help='Root directory to process')
    parser.add_argument('--dry-run', action='store_true',
                        help='Show what would change without modifying files')
    parser.add_argument('--ext', nargs='+', metavar='EXT',
                        help='Only process files with these extensions (e.g. .txt .md)')
    args = parser.parse_args()

    root = Path(args.directory)
    if not root.is_dir():
        print(f"Error: '{root}' is not a directory.", file=sys.stderr)
        sys.exit(1)

    extensions = [e if e.startswith('.') else f'.{e}' for e in args.ext] if args.ext else []

    print(f"Processing: {root.resolve()}")
    if args.dry_run:
        print("Mode: DRY RUN (no files will be changed)")
    if extensions:
        print(f"Extensions: {', '.join(extensions)}")
    print()

    process_directory(root, dry_run=args.dry_run, extensions=extensions)


if __name__ == '__main__':
    main()
