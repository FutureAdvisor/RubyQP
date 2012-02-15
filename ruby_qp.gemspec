# -*- encoding: utf-8 -*-
$:.push File.expand_path("../lib", __FILE__)
require "ruby_qp/version"

Gem::Specification.new do |s|
  s.name        = "ruby_qp"
  s.version     = RubyQp::VERSION
  s.authors     = ["Josh Tokle"]
  s.email       = ["josh@futureadvisor.com"]
  s.homepage    = ""
  s.summary     = %q{CQP wrapper.}
  s.description = %q{CQP wrapper.}

  s.rubyforge_project = "ruby_qp"

  s.files         = `git ls-files`.split("\n")
  s.test_files    = `git ls-files -- {test,spec,features}/*`.split("\n")
  s.extensions    << 'ext/extconf.rb'
  s.executables   = `git ls-files -- bin/*`.split("\n").map{ |f| File.basename(f) }
  s.require_paths = ["lib"]
end
