# Copyright (c) 2021 Ed Harry, Wellcome Sanger Institute, Genome Research Limited
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

project(
  'LinkStats',
  'cpp',
default_options : ['cpp_std=c++20'],
license : 'MIT',
meson_version : '>=0.60.0',
version : '0.0.1-1-ge17d43d'
)


name = meson.project_name()

lib_src = name + '.cpp'
lib_header = name + '.hpp'

test_src = 'test.cpp'
test_file = 'test.bam'

py_CExt_name = '_LinkStats_C'
py_CExt_src = 'LinkStats_PyCExt.cpp'


flags = ['-Ofast']
if get_option('buildtype') == 'debug'
    flags = ['-O0', '-g', '-DDEBUG']
endif


deps = [  dependency('threads'),
          dependency('htslib', version : '>= 1.14')]
cpp_dep = meson.get_compiler('cpp').find_library('stdc++', required : true)


main_lib = static_library(name, sources : lib_src,
                                dependencies : deps + cpp_dep,
                                cpp_args : flags)

test(name + '_test',  executable(name,  test_src,
                                        dependencies : cpp_dep, 
                                        link_with : main_lib, 
                                        cpp_args : flags),
                      args : meson.project_source_root() + '/' + test_file)


py_inst = import('python').find_installation('python3')

py_inst.extension_module(py_CExt_name,  py_CExt_src, 
                                        dependencies : [py_inst.dependency()] + cpp_dep,
                                        link_with : main_lib,
                                        cpp_args : flags,
                                        install : true)
