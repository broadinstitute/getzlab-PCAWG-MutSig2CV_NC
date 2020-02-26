function Y = move_field_to_after(X,fld1,fld2)
% Y = move_field_to_after(X,fld1,fld2)

demand_fields(X,{fld1,fld2});
f = fieldnames(X);

for i=1:length(f)
  if strcmp(f{i},fld1), continue; end
  Y.(f{i}) = X.(f{i});
  if strcmp(f{i},fld2)
    Y.(fld1) = X.(fld1);
  end
end
