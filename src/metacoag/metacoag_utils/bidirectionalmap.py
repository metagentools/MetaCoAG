#!/usr/bin/env python3

"""An invertible dictionary implementation.

Adapted from https://stackoverflow.com/questions/3318625/.
"""


class BidirectionalError(Exception):
    """Must set a unique value in a BijectiveMap."""

    def __init__(self, value):
        self.value = value
        message = f'The value "{value}" is already in the mapping.'
        super().__init__(message)


class BidirectionalMap(dict):
    """Invertible map."""

    def __init__(self, inverse=None):
        if inverse is None:
            inverse = self.__class__(inverse=self)
        self.inverse = inverse

    def __setitem__(self, key, value):
        if value in self.inverse:
            raise BidirectionalError(value)

        self.inverse._set_item(value, key)
        self._set_item(key, value)

    def __delitem__(self, key):
        self.inverse._del_item(self[key])
        self._del_item(key)

    def _del_item(self, key):
        super().__delitem__(key)

    def _set_item(self, key, value):
        super().__setitem__(key, value)
