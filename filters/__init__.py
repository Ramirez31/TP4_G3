from importlib import import_module

from .base_filter import base_filter

def create(filter_name, *args, **kwargs):

    try:
        if '.' in filter_name:
            module_name, class_name = filter_name.rsplit('.', 1)
        else:
            module_name = filter_name
            class_name = filter_name.capitalize()
        filter_module = import_module('.' + module_name, package='filters') 
        filter_class = getattr(filter_module, class_name)
        instance = filter_class(*args, **kwargs)

    except (AttributeError, ModuleNotFoundError):
        raise ImportError('{} is not part of our filter selection!'.format(filter_name))
    else:
        if not issubclass(filter_class, base_filter):
            raise ImportError("We currently don't have {}, but you are welcome to send in the request for it!".format(filter_class))

    return instance