"""Tests for kmtools.utils."""

import time
from io import StringIO
from unittest.mock import patch

from kmtools.utils import Utils


class TestLog:
    def test_log_verbose_prints_to_stderr(self):
        with patch("sys.stderr", new_callable=StringIO) as mock_stderr:
            Utils.log("test message", verbose=True)
            assert "test message" in mock_stderr.getvalue()

    def test_log_not_verbose_prints_nothing(self):
        with patch("sys.stderr", new_callable=StringIO) as mock_stderr:
            Utils.log("test message", verbose=False)
            assert mock_stderr.getvalue() == ""

    def test_log_default_verbose_is_false(self):
        with patch("sys.stderr", new_callable=StringIO) as mock_stderr:
            Utils.log("test message")
            assert mock_stderr.getvalue() == ""


class TestTimeIt:
    def test_time_it_returns_function_result(self):
        result = Utils.time_it("test", lambda x: x * 2, 5, verbose=False)
        assert result == 10

    def test_time_it_verbose_logs_timing(self):
        with patch("sys.stderr", new_callable=StringIO) as mock_stderr:
            Utils.time_it("my_label", lambda: 42, verbose=True)
            output = mock_stderr.getvalue()
            assert "my_label" in output
            assert "completed in" in output

    def test_time_it_passes_kwargs(self):
        def func(a, b=10):
            return a + b

        result = Utils.time_it("test", func, 5, b=20, verbose=False)
        assert result == 25
