require 'rubygems'
require 'rubygems/package_task'
require 'rake/clean'

EXT_CONF = 'ext/extconf.rb'
MAKEFILE = 'ext/Makefile'
MODULE = 'ext/ruby_qp.so'
SRC = Dir.glob('ext/*.c')
SRC << MAKEFILE

CLEAN.include [ 'ext/*.o', MODULE, 'pkg' ]
CLOBBER.include [ 'config.save', 'ext/mkmf.log', 'ext/extconf.h',
                  MAKEFILE ]

PKG_FILES = FileList[
    "Rakefile", "ext/*.[ch]", "ext/extconf.rb",
]

spec = Gem::Specification.new do |s|
  s.platform = Gem::Platform::RUBY
  s.summary = "Ruby wrapper for the cqp library."
  s.name = "ruby_qp"
  s.version = "0.0.1"
  s.requirements << 'GSL (GNU Scientific Library)'
  s.require_paths << 'ext'
  s.extensions = "ext/extconf.rb"
  s.files = PKG_FILES
  s.description = "ruby wrapper for the cqp library"
  s.author = "Joshua Tokle"
  s.email = "josh@futureadvisor.com"
  s.homepage = "http://www.futureadvisor.com"
end

file MAKEFILE => EXT_CONF do |t|
    Dir::chdir(File::dirname(EXT_CONF)) do
        unless sh "ruby #{File::basename(EXT_CONF)}"
            $stderr.puts "Failed to run extconf"
            break
        end
    end
end
file MODULE => SRC do |t|
    Dir::chdir(File::dirname(EXT_CONF)) do
        unless sh "make"
            $stderr.puts "make failed"
            break
        end
    end
end
desc "Build the native library"
task :build => MODULE

Gem::PackageTask.new(spec) do |pkg|
    pkg.need_tar = true
    pkg.need_zip = true
end
