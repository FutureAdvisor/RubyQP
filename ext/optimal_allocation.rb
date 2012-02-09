require 'rubyqp'

def optimal_allocation(funds, target_allocation, sector_weights = Array.new(target_allocation.size, 1))
  fund_ids = funds.keys
  sector_ids = target_allocation.keys

  # In order to work with the optimization method, put sectors in rows and funds in columns
  fund_matrix = []
  sector_ids.each do |sector|
    row = []
    fund_ids.each do |fund|
      row << funds[fund][sector]
    end
    fund_matrix << row
  end 

  fund_matrix
end
