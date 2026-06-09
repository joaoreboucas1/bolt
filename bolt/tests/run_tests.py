"""
    Bolt test suite
    Usage: `make test` or `python3 ./tests/run_tests.py`
"""

import subprocess
import sys
from time import perf_counter
from enum import Enum, auto

class TestResult(Enum):
    FAIL_COMPILE_C = auto()
    FAIL_RUN_C = auto()
    UNEXPECTED_OUTPUT_C = auto()
    FAIL_RUN_PY = auto()
    UNEXPECTED_OUTPUT_PY = auto()
    SUCCESS = auto()

class Test():
    def __init__(self, name, expected_output=None):
        self.name = name
        self.expected_output = expected_output

    def run(self):
        proc = subprocess.run([
            "gcc",
            "-o", f"tests/{self.name}",
            f"tests/{self.name}.c",
            "build/dverk.o",
            "-L/opt/homebrew/lib/gcc/current",
            "-lm", "-lgsl", "-lgslcblas", "-ggdb", "-lgfortran"
        ])
        if proc.returncode != 0:
            return TestResult.FAIL_COMPILE_C, None
        
        start = perf_counter()
        proc = subprocess.run([f"./tests/{self.name}"], capture_output=True)
        elapsed_time_c = perf_counter() - start
        if proc.returncode != 0:
            return TestResult.FAIL_RUN_C, None
        
        if self.expected_output is not None and proc.stdout.decode("utf-8") != self.expected_output:
            return TestResult.UNEXPECTED_OUTPUT_C, None
        
        start = perf_counter()
        proc = subprocess.run(["python3", f"./tests/{self.name}.py"], capture_output=True)
        elapsed_time_py = perf_counter() - start
        if proc.returncode != 0:
            return TestResult.FAIL_RUN_PY, None

        if self.expected_output is not None and proc.stdout.decode("utf-8") != self.expected_output:
            return TestResult.UNEXPECTED_OUTPUT_PY, None

        return TestResult.SUCCESS, (elapsed_time_c, elapsed_time_py)    

tests = [
    Test("distances", expected_output="z = 0.50, \\chi = 1960.31 Mpc || z = 1.00, \\chi = 3413.50 Mpc\n"),
    Test("thermo", expected_output="At z = 1300, \\kappa'(z) = 0.16246, \\kappa(z) = 1.44999, g(z)/H_0 = 170.28\n"),
    Test("matter_tk")
]

if __name__ == "__main__":
    for test in tests:
        result, times = test.run()
        if result == TestResult.SUCCESS:
            elapsed_time_c, elapsed_time_py = times
            print(f"Test {test.name:12}=> {result} (time_c = {elapsed_time_c:.8f}s, time_py = {elapsed_time_py:.8f}s)")
        else: 
            print(f"Test {test.name:12}=> {result}")