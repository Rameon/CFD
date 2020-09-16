%% CALCULATE RESIDUAL OF GRID %%
function resid = calcResid(V, V_free, normals, areas, nci, ncj)

nE = normals{1};
nN = normals{2};
nW = normals{3};
nS = normals{4};

sE = areas{1};
sN = areas{2};
sW = areas{3};
sS = areas{4};

e_flux = cell(nci,ncj);
w_flux = cell(nci,ncj);
n_flux = cell(nci,ncj);
s_flux = cell(nci,ncj);

resid = cell(nci,ncj);


free_vel = [V_free(2) V_free(3)];
free_a   = speedsound(V_free(4), V_free(1));

for i = 1:nci
    for j = 1:ncj
        cell_a = speedsound(V{i,j}(4),V{i,j}(1));
        vel    = [V{i,j}(2) V{i,j}(3)];
        conM = dot(vel,nE{i,j})/cell_a;
        
        % East face is at the outflow boundary
        if i == nci
            nB = nE{i,j};
            if conM >= 1
                Vb = V{i,j};
                e_flux{i,j} = flux(Vb,nB);
            else
                % Calculate positive, negative Riemann invariants
                rpos = (vel.*nB) + 2/(1.4-1)*cell_a;
                rneg = (free_vel.*nB) - 2/(1.4-1)*free_a;
                un = (rpos+rneg)/2;
                ub = vel + (un - vel).*nB;                 

                Vb = [V{i,j}(1) ub(1) ub(2) V{i,j}(4)];
                e_flux{i,j} = flux(Vb,nB);  
            end

        else
            if conM >= 1
                e_flux{i,j} = flux(V{i,j},nE{i,j});
            else
                e_flux{i,j} = fpos(V{i,j},nE{i,j}) + fneg(V{i+1,j},nE{i,j});
            end
        end
        
        
        conM = dot(vel,nN{i,j})/cell_a;
        % North face is invicid wall
        if j == ncj
            nB = nN{i,j};
            conV = dot(vel, nB);
            velB = vel - conV*nB;
            
            Vb = [V{i,j}(1) velB(1) velB(2) V{i,j}(4)];
            n_flux{i,j} = flux(Vb,nB);
        else
            if conM >=1
                n_flux{i,j} = flux(V{i,j},nN{i,j});
            else
                n_flux{i,j} = fpos(V{i,j},nN{i,j}) + fneg(V{i,j+1},nN{i,j});
            end
        end
        
        
        conM = dot(vel,nW{i,j})/cell_a;
        % West face is inflow boundary
        if i == 1
            nB = nW{i,j}; 
            if conM >= 1
                Vb = V_free;
                w_flux{i,j} = flux(Vb,nB);
            else
                % Calculate positive, negative Riemann invariant
                rpos = (vel.*nB) + 2/(1.4-1)*cell_a;
                rneg = (free_vel.*nB) - 2/(1.4-1)*free_a;
                un = (rpos+rneg)/2;
                ub = free_vel + (un + free_vel).*nB;
                
                Vb = [V_free(1) ub(1) ub(2) V_free(4)];

                w_flux{i,j} = flux(Vb,nB);
            end
        else
            if conM >=1
                w_flux{i,j} = flux(V{i,j},nW{i,j});
            else
                w_flux{i,j} = -1*e_flux{i-1,j};
            end
        end
        
        
        conM = dot(vel,nS{i,j})/cell_a;
        
        % South face is invicid wall
        if j == 1
            nB = nS{i,j};
            conV = dot(vel, nB);
            velB = vel - conV*nB;
            
            Vb = [V{i,j}(1) velB(1) velB(2) V{i,j}(4)];
            
            s_flux{i,j} = flux(Vb,nB);
        else
            if conM >= 1
                s_flux{i,j} = flux(V{i,j},nS{i,j});
            else
                s_flux{i,j} = -1*n_flux{i,j-1};                            
            end
        end
        

        resid{i,j} = e_flux{i,j}*sE(i,j) + n_flux{i,j}*sN(i,j) + w_flux{i,j}*sW(i,j) + s_flux{i,j}*sS(i,j);
    end
end


