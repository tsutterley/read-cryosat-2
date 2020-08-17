import pytest

def pytest_addoption(parser):
    parser.addoption("--username", action="store", help="ESA FTP username")
    parser.addoption("--password", action="store", help="ESA FTP password")

@pytest.fixture
def username(request):
    """ Returns ESA FTP username """
    return request.config.getoption("--username")

@pytest.fixture
def password(request):
    """ Returns ESA FTP password """
    return request.config.getoption("--password")
