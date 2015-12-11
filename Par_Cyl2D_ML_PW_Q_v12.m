function Par_Cyl2D_ML_PW_Q_v12(filename,N_cur)
%code to evaluate integral from derivation of spiral field value as a
%function position
if ischar(N_cur)
    N_cur=str2num(N_cur);
end

load([filename '.mat'])

N_list=N_dist(N_cur).N;

%try and catch block to allow it to operate with different number of layers
try
    try
        for rr1=1:length(r1)
            for rr2=1:length(r2)
                for rr3=1:length(r3)
                    for jj=1:length(N_list)
                        [Q(jj,rr1,rr2,rr3,:)]=Cyl2D_ML_PW_Q_v12(N_list(jj),cumsum([r1(rr1) r2(rr2) r3(rr3)]),materials,isEz);
                    end
                end
            end
        end
    catch
        for rr1=1:length(r1)
            for rr2=1:length(r2)
                for jj=1:length(N_list)
                    [Q(jj,rr1,rr2,:)]=Cyl2D_ML_PW_Q_v12(N_list(jj),cumsum([r1(rr1) r2(rr2)]),materials,isEz);
                end
            end
        end
    end
catch
    for rr1=1:length(r1)
        for jj=1:length(N_list)
            [Q(jj,rr1,:)]=Cyl2D_ML_PW_Q_v12(N_list(jj),cumsum([r1(rr1)]),materials,isEz);
        end
    end
end

save([filename '_out' num2str(N_cur) '.mat'],'Q')
