"""
    Bolt test suite
    Usage: `make test` or `python3 ./tests/run_tests.py`
"""

import subprocess
from enum import Enum, auto

class TestResult(Enum):
    FAIL_COMPILE_C = auto()
    FAIL_RUN_C = auto()
    UNEXPECTED_OUTPUT_C = auto()
    FAIL_RUN_PY = auto()
    UNEXPECTED_OUTPUT_PY = auto()
    SUCCESS = auto()

def test_background():
    proc = subprocess.run([
        "gcc",
        "-o", "tests/distances",
        "tests/distances.c",
        "-lm", "-lgsl", "-lgslcblas", "-ggdb"
    ])
    if proc.returncode != 0:
        return TestResult.FAIL_COMPILE_C
    
    expected = "z = 0.50, \\chi = 1960.31 Mpc || z = 1.00, \\chi = 3413.50 Mpc\n"
    proc = subprocess.run(["./tests/distances"], capture_output=True)
    if proc.returncode != 0:
        return TestResult.FAIL_RUN_C
    
    if not proc.stdout.decode("utf-8") == expected:
        return TestResult.UNEXPECTED_OUTPUT_C
    
    proc = subprocess.run(["python3", "./tests/distances.py"], capture_output=True)
    if proc.returncode != 0:
        return TestResult.FAIL_RUN_PY

    if not proc.stdout.decode("utf-8") == expected:
        return TestResult.UNEXPECTED_OUTPUT_PY

    return TestResult.SUCCESS

def test_matter_tk():
    proc = subprocess.run([
        "gcc",
        "-o", "tests/matter_tk",
        "tests/matter_tk.c",
        "-lm", "-lgsl", "-lgslcblas", "-ggdb"
    ])

    if proc.returncode != 0: return TestResult.FAIL_COMPILE_C

    proc = subprocess.run(["./tests/matter_tk"], capture_output=True)
    if proc.returncode != 0:
        return TestResult.FAIL_RUN_C
    
    proc = subprocess.run(["python3", "./tests/matter_tk.py"], capture_output=True)
    if proc.returncode != 0:
        return TestResult.FAIL_RUN_PY

    return TestResult.SUCCESS

tests = {
    "comoving distances": test_background,
    "matter transfer function": test_matter_tk
}

if __name__ == "__main__":
    for name, test in tests.items():
        print(f"Test {name}: {test()}")