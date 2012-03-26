module RubyQp
  class OptimizationError < StandardError; end

  class InfeasibleProblemDetectedError < OptimizationError; end
  class SearchDirectionBecomesTooSmallError < OptimizationError; end
  class DivergingIteratesError < OptimizationError; end
  class MaximumIterationsExceededError < OptimizationError; end
  class RestorationFailedError < OptimizationError; end
  class ErrorInStepComputationError < OptimizationError; end
  class InvalidOptionError < OptimizationError; end
  class NotEnoughDegreesOfFreedomError < OptimizationError; end
  class InvalidProblemDefinitionError < OptimizationError; end
  class UnrecoverableExceptionError < OptimizationError; end
  class NonIpoptExceptionThrownError < OptimizationError; end
  class InternalError < OptimizationError; end
end
