function [first_item last_item] = calculate_work_splits(num_splits, num_items)

if num_splits > num_items
  nextra = num_splits - num_items;
  num_splits = num_items;
end

tmp = round([1:num_items/num_splits:num_items]);
tmp(num_splits+1) = num_items+1;

first_item = tmp(1:num_splits);
last_item = tmp(2:num_splits+1)-1;

if exist('nextra','var')
  first_item = [first_item zeros(1,nextra)];
  last_item = [last_item -ones(1,nextra)];
end

