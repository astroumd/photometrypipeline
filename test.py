#!/usr/bin/env python

import os

from zipfile import ZipFile

from reduction import reduction


def test(noplot=True):
    print('starting test')
    base_path = os.path.abspath(os.path.dirname(__file__))  # /path/to/photometrypipeline
    print('base path', base_path)
    test_path = os.path.join(base_path, 'test')  # /path/to/photometrypipeline/test
    print('test_path', test_path)
    input_data_path = os.path.join(test_path, 'extracted_data')  # /path/to/photometrypipeline/test/extracted_data
    zip_file = os.path.join(test_path, 'test.zip')  # /path/to/photometrypipeline/test/test.zip
    print('zip_file', zip_file)
    zf = ZipFile(zip_file, 'r')
    print('extracting files')
    zf.extractall(input_data_path)
    print('extraction complete')
    reduction(input_data_path, noplot=True)


if __name__ == '__main__':
    test()
