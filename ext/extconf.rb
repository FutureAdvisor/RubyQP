require 'mkmf'

extension_name = 'ruby_qp'

find_executable('pkg-config') or raise 'pkg-config should be installed'
gsl_vars = pkg_config('gsl') or raise 'GSL not found!'
ipopt_vars = pkg_config('ipopt') or raise 'Ipopt not found!'

create_header
create_makefile(extension_name)
