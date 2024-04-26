function pca_dim(explained,btn_r,btn_nr,btn_w,btn_noise)
    b = bar(explained(1:20));
    for i =1:length(btn_r)
        if btn_r(i).Value
            b.CData(i,:) = [.8 0 .2];
        end
    end
    for i =1:length(btn_nr)
        if btn_nr(i).Value
            b.CData(i,:) = [250/255 225/255 90/255];
        end
    end
    for i =1:length(btn_w)
        if btn_w(i).Value
            b.CData(i,:) = [0 0.2 .8];
        end
    end
    for i =1:length(btn_noise)
        if btn_noise(i).Value
            b.CData(i,:) = [1 1 1];

        end
    end

end