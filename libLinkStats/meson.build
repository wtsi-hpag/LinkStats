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
main_src = 'main.cpp'

flags = ['-Ofast']
if get_option('buildtype') == 'debug'
    flags = ['-O0', '-g', '-DDEBUG']
endif
flags += ['-DPV=' + meson.project_version()]

deps = [dependency('threads'), dependency('htslib', version : '>= 1.14')]
cpp_dep = meson.get_compiler('cpp').find_library('stdc++', required : true)

main_lib = static_library(name, sources : lib_src, dependencies : deps + cpp_dep, cpp_args : flags, install: true)
install_headers(lib_header)
executable(name, main_src, dependencies : cpp_dep, link_with : main_lib, cpp_args : flags)
