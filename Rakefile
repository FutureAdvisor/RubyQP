require 'rake/clean'
require 'rake/gempackagetask'

EXT_CONF = 'ext/extconf.rb'
MAKEFILE = 'ext/Makefile'
MODULE = 'ext/rubyqp.so'
SRC = Dir.glob('ext/*.c')
SRC << MAKEFILE

CLEAN.include [ 'ext/*.o', MODULE, 'pkg' ]
CLOBBER.include [ 'config.save', 'ext/mkmf.log', 'ext/extconf.h',
                  MAKEFILE ]

PKG_FILES = FileList[
    "Rakefile", "ext/*.[ch]", "ext/extconf.rb",
]

SPEC = Gem::Specification.new do |s|
    s.name = "rubyqp"
    s.version = "0.0.1"
    s.email = "josh@futureadvisor.com"
    s.homepage = "http://futureadvisor.com/"
    s.summary = "ruby wrapper for the cqp library"
    s.files = PKG_FILES
    s.required_ruby_version = '>= 1.8.7'
    s.extensions = "ext/extconf.rb"
    s.author = "Joshua Tokle"
    s.rubyforge_project = "None"
    s.description = "ruby wrapper for the cqp library"
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

Rake::GemPackageTask.new(SPEC) do |pkg|
    pkg.need_tar = true
    pkg.need_zip = true
end
