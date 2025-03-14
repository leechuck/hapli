#!/usr/bin/env python3
"""
This script is a wrapper around hapli.tools.gfa_to_gff3

This file is kept for backward compatibility but will be removed in a future version.
Please use the module directly: python -m hapli.tools.gfa_to_gff3
"""

import os
import sys
import warnings

# Add parent directory to path to allow importing from hapli
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '.')))

# Show deprecation warning
warnings.warn(
    "This script is now a wrapper around hapli.tools.gfa_to_gff3. "
    "Please use the module directly: python -m hapli.tools.gfa_to_gff3",
    DeprecationWarning,
    stacklevel=2
)

# Import and run the main function from the module
from hapli.tools.gfa_to_gff3 import main

if __name__ == "__main__":
    sys.exit(main())
