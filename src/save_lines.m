function save_lines(L,fname)

L = as_column(L);

X{1} = L;
X{2} = repmat({char(10)},length(L),1);
X = strcat(X{:});
X = [X{:}];
save_textfile(X,fname);
