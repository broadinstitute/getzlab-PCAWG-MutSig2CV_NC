function [first_item last_item] = calculate_work_split(split_num, num_splits, num_items)
% [first_item last_item] = calculate_work_split(split_num, num_splits, num_items)
%
% Mike Lawrence
% 2008-06-13

if split_num<1 || split_num>num_splits
  first_item = 0;
  last_item = -1;
else
  [f l] = calculate_work_splits(num_splits,num_items);
  first_item = f(split_num);
  last_item = l(split_num);
end
