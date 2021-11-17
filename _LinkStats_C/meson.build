project(
  'LinkStats',
  'cpp',
default_options : ['cpp_std=c++20'],
license : 'MIT',
meson_version : '>=0.60.0',
version : '0.0.1'
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
