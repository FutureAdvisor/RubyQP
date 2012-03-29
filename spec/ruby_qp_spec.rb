require 'rubygems'
require 'ruby_qp'

describe RubyQp do

  # Parameters to the quadratic programming problem
  #   ||Ax - b||
  # under constraints
  #   x_lower <= x <= x_upper
  #   g_lower <= Gx <= g_upper
  before do

    @a_mat = [[0.5315652, 0, 1], 
              [0, 1, 0], 
              [0, 0, 0], 
              [0.1287067, 0, 0]]
    @b_vec = [0.14742857, 0.07371429, 0.14742857, 0.07371429]

    @x_init = [0.33, 0.33, 0.34]
    @x_lower = [0, 0, 0]
    @x_upper = [1, 1, 1]

    @g_mat = [[1, 1, 1]]
    @g_lower = [1]
    @g_upper = [1]

    # solution to the QP problem as computed in R
    @x_vec = [0.7806328, 0.2193672, 0.0]
  end


  describe '#solve_dist_full' do
    it 'returns a Hash' do
      h = RubyQp::solve_dist_full @a_mat, @b_vec, @x_lower, @x_upper, @g_mat, @g_lower, @g_upper, @x_init
      h.should be_an_instance_of(Hash)
    end

    it 'solves quadratic programming problems' do
      h = RubyQp::solve_dist_full @a_mat, @b_vec, @x_lower, @x_upper, @g_mat, @g_lower, @g_upper, @x_init
      x_vec = h[:solution]
      x_vec.each_index do |ix|
        x_vec[ix].should be_close(@x_vec[ix], 1e-7)
      end
    end
  end

  describe '#solve_dist' do
    it 'returns an Array' do
      ary = RubyQp::solve_dist @a_mat, @b_vec, @x_lower, @x_upper, @g_mat, @g_lower, @g_upper, @x_init
      ary.should be_an_instance_of(Array)
    end

    it 'solves min-distance problems' do
      x_vec = RubyQp::solve_dist @a_mat, @b_vec, @x_lower, @x_upper, @g_mat, @g_lower, @g_upper, @x_init
      x_vec.each_index do |ix|
        x_vec[ix].should be_close(@x_vec[ix], 1e-7)
      end
    end
  end

end
