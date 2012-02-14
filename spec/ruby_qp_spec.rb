require 'rubygems'
require 'ruby_qp'

describe RubyQp do

  # Parameters to the quadratic programming problems
  #   (1/2)(x^t)Qx + (q^t)x
  # and
  #   ||Mx - m||
  # under constraints
  #   Ax = b
  #   Cx >= d
  # Parameters Q, q, M, and m are chosen so that both problems are solved by
  # @x_vec.
  before do
    @q_mat = [[0.2991269, 0, 0.5315652],
              [0,         1, 0],
              [0.5315652, 0, 1]]
    @q_vec = [-0.08785541, -0.07371429, -0.14742857]
    @a_mat = [[1, 1, 1]]
    @b_vec = [1]
    @c_mat = [[1, 0, 0],
              [0, 1, 0],
              [0, 0, 1]]
    @d_vec = [0, 0, 0]

    @m_mat = [[0.5315652, 0, 1], 
              [0, 1, 0], 
              [0, 0, 0], 
              [0.1287067, 0, 0]]
    @m_vec = [0.14742857, 0.07371429, 0.14742857, 0.07371429]

    # solution to the QP problem as computed in R
    @x_vec = [0.7806328, 0.2193672, 0.0]
    @lagrange_eq = [0.1456529]
    @lagrange_ineq = [0.0, 0.0, 0.1218758]
  end


  describe '#solve_full' do
    it 'returns a Hash' do
      h = RubyQp::solve_full @q_mat, @q_vec, @a_mat, @b_vec, @c_mat, @d_vec
      h.should be_an_instance_of(Hash)
    end

    it 'solves quadratic programming problems' do
      h = RubyQp::solve_full @q_mat, @q_vec, @a_mat, @b_vec, @c_mat, @d_vec
      x_vec = h["solution"]
      x_vec.each_index do |ix|
        x_vec[ix].should be_close(@x_vec[ix], 1e-7)
      end
    end

    it 'computes the Lagrangian of its constraints' do 
      h = RubyQp::solve_full @q_mat, @q_vec, @a_mat, @b_vec, @c_mat, @d_vec

      lagrange_eq = h["lagrange_eq"]
      lagrange_eq.each_index do |ix|
        lagrange_eq[ix].should be_close(@lagrange_eq[ix], 1e-7)
      end

      lagrange_ineq = h["lagrange_ineq"]
      lagrange_ineq.each_index do |ix|
        lagrange_ineq[ix].should be_close(@lagrange_ineq[ix], 1e-7)
      end
    end
  end

  describe '#solve' do
    it 'returns an Array' do
      ary = RubyQp::solve @q_mat, @q_vec, @a_mat, @b_vec, @c_mat, @d_vec
      ary.should be_an_instance_of(Array)
    end

    it 'solves quadtratic programming problems' do
      x_vec = RubyQp::solve @q_mat, @q_vec, @a_mat, @b_vec, @c_mat, @d_vec
      x_vec.each_index do |ix|
        x_vec[ix].should be_close(@x_vec[ix], 1e-7)
      end
    end
  end

  describe '#solve_dist_full' do
    it 'returns a Hash' do
      h = RubyQp::solve_dist_full @m_mat, @m_vec, @a_mat, @b_vec, @c_mat, @d_vec
      h.should be_an_instance_of(Hash)
    end

    it 'solves quadratic programming problems' do
      h = RubyQp::solve_dist_full @m_mat, @m_vec, @a_mat, @b_vec, @c_mat, @d_vec
      x_vec = h["solution"]
      x_vec.each_index do |ix|
        x_vec[ix].should be_close(@x_vec[ix], 1e-7)
      end
    end

    it 'computes the Lagrangian of its constraints' do 
      h = RubyQp::solve_dist_full @m_mat, @m_vec, @a_mat, @b_vec, @c_mat, @d_vec

      lagrange_eq = h["lagrange_eq"]
      lagrange_eq.each_index do |ix|
        lagrange_eq[ix].should be_close(@lagrange_eq[ix], 1e-7)
      end

      lagrange_ineq = h["lagrange_ineq"]
      lagrange_ineq.each_index do |ix|
        lagrange_ineq[ix].should be_close(@lagrange_ineq[ix], 1e-7)
      end
    end
  end

  describe '#solve_dist' do
    it 'returns an Array' do
      ary = RubyQp::solve_dist @m_mat, @m_vec, @a_mat, @b_vec, @c_mat, @d_vec
      ary.should be_an_instance_of(Array)
    end

    it 'solves min-distance problems' do
      x_vec = RubyQp::solve_dist @m_mat, @m_vec, @a_mat, @b_vec, @c_mat, @d_vec
      x_vec.each_index do |ix|
        x_vec[ix].should be_close(@x_vec[ix], 1e-7)
      end
    end
  end

end
