#!/usr/bin/env python3
"""
This script has been moved to hapli.tools.generate_test_data

This file is kept for backward compatibility but will be removed in a future version.
Please update your scripts to use the new location.
"""

import os
import sys
import warnings

# Add parent directory to path to allow importing from hapli
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

# Show deprecation warning
warnings.warn(
    "This script has been moved to hapli.tools.generate_test_data. "
    "Please update your scripts to use the new location.",
    DeprecationWarning,
    stacklevel=2
)

# Import and run the main function from the new location
from hapli.tools.generate_test_data import main

if __name__ == "__main__":
    sys.exit(main())
