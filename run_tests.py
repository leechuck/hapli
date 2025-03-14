#!/usr/bin/env python3
"""
Test runner for hapli.

This script runs all tests for the hapli project.
"""

import os
import sys
import unittest
import argparse
import logging

def setup_logging(debug=False):
    """Configure logging based on debug flag."""
    log_level = logging.DEBUG if debug else logging.INFO
    logging.basicConfig(
        level=log_level,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )

def discover_and_run_tests(test_dir=None, pattern='test_*.py', verbosity=1):
    """Discover and run tests in the specified directory."""
    if test_dir is None:
        # Use the directory where this script is located
        test_dir = os.path.dirname(os.path.abspath(__file__))
        # If run from project root, use the tests directory
        if os.path.basename(test_dir) != 'tests':
            test_dir = os.path.join(test_dir, 'tests')
    
    logging.info(f"Discovering tests in {test_dir} with pattern {pattern}")
    
    # Discover tests
    loader = unittest.TestLoader()
    suite = loader.discover(test_dir, pattern=pattern)
    
    # Run tests
    runner = unittest.TextTestRunner(verbosity=verbosity)
    result = runner.run(suite)
    
    return result

def run_specific_test(test_path, verbosity=1):
    """Run a specific test file."""
    logging.info(f"Running test file: {test_path}")
    
    # Add the directory to sys.path if needed
    test_dir = os.path.dirname(os.path.abspath(test_path))
    if test_dir not in sys.path:
        sys.path.insert(0, test_dir)
    
    # Load the test module
    test_module = os.path.basename(test_path)
    if test_module.endswith('.py'):
        test_module = test_module[:-3]
    
    # Import the module and run tests
    try:
        module = __import__(test_module)
        suite = unittest.TestLoader().loadTestsFromModule(module)
        runner = unittest.TextTestRunner(verbosity=verbosity)
        result = runner.run(suite)
        return result
    except ImportError as e:
        logging.error(f"Error importing test module {test_module}: {e}")
        return None

def main():
    """Main function to run tests."""
    parser = argparse.ArgumentParser(description='Run tests for hapli')
    parser.add_argument('--test-dir', help='Directory containing tests')
    parser.add_argument('--pattern', default='test_*.py', help='Pattern to match test files')
    parser.add_argument('--test-file', help='Run a specific test file')
    parser.add_argument('--verbosity', type=int, default=2, help='Verbosity level (1-3)')
    parser.add_argument('--debug', action='store_true', help='Enable debug output')
    
    args = parser.parse_args()
    
    # Setup logging
    setup_logging(debug=args.debug)
    
    # Run tests
    if args.test_file:
        result = run_specific_test(args.test_file, verbosity=args.verbosity)
    else:
        result = discover_and_run_tests(args.test_dir, args.pattern, verbosity=args.verbosity)
    
    # Return appropriate exit code
    if result and result.wasSuccessful():
        logging.info("All tests passed!")
        return 0
    else:
        logging.error("Some tests failed.")
        return 1

if __name__ == '__main__':
    sys.exit(main())
