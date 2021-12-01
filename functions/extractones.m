function out = extractones(ins)

out = {};
insnum = find(ins);
diffnum = (diff([insnum(1);insnum]) > 1);
st = insnum(1);
k = 1;
for i = 1:numel(diffnum)
    if diffnum(i) > 0
        out{k} = st:insnum(i-1);
        st = insnum(i);
        k = k + 1;
    end
end
         

end