from common import fill_values


def test_fill_func():
    assert fill_values([None, 1, 2]) == [2, 1, 2]
    assert fill_values([10, None, 5]) == [10, 2, 5]
    assert fill_values([10, 2, None]) == [10, 2, 5]
