"""Pytest configuration for the leakdecode suite."""


def pytest_configure(config):
    config.addinivalue_line(
        "markers",
        "slow: integration test that runs the quick benchmark (still collected by default)",
    )
