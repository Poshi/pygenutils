#!/bin/bash

. .venv/bin/activate
python3 -m build
python3 -m pip install --upgrade "$(find dist -name "*.whl" -printf "%T@ %p\n" | sort -rn | head -1 | cut -f2 -d' ')"