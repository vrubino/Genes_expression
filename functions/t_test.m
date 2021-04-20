function [p, pFDR, h] = t_test(E, mask, a)
p = zeros(1,size(mask,2));
    for k = 1:size(mask,2)
        tmp_mask = mask(:,k);
        [~,p(k)] = ttest2(E(:,tmp_mask), E(:,~tmp_mask), ...
            'tail', 'right', 'Vartype','unequal');
    end
pFDR = mafdr(p);
h = pFDR < a;