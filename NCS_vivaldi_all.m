% F. Dabek, R. Cox, and F. Kaashoek. Vivaldi: A Decentralized Network Coordinate System. In Proc. of ACM SIGCOMM, 2004.

function [predicted_matrix] = NCS_vivaldi_all(DATA, dim, N, K, option)
% option = 0, basic Vivaldi (Dabek et al., SIGCOMM'04)
% option = 1, basic Vivaldi with height element (Dabek et al., SIGCOMM'04)
% option = 2, TIV aware Vivaldi (Wang et al., IMC'07)
    if (option == 0)
        [coord_all] = NCS_vivaldi_basic(DATA, dim, length(DATA), K);   
        predicted_matrix = (zeros(length(DATA), length(DATA)));
        for x = 1:length(DATA)        
            for y = 1:length(DATA)
                predicted_matrix(x,y) = Dist(coord_all(x,:), coord_all(y,:)); 
            end
        end    
    elseif (option == 1);
        [coord_all] = NCS_vivaldi_h(DATA, dim, length(DATA), K);    
%        total_vivaldi_fperror = [total_vivaldi_fperror vivaldi_fperror];
        for x = 1:length(DATA)
            %coord_all(x, :)
            for y = 1:length(DATA)
                predicted_matrix(x,y) = Dist_h(coord_all(x,:), coord_all(y,:)); 
            end;
        end;    
        
    elseif (option == 2);
        [coord_all, vivaldi_fperror] = NCS_vivaldi_tivaware(DATA, dim, length(DATA), K, 0);   
        predicted_matrix = (zeros(length(DATA), length(DATA)));
        for x = 1:length(DATA)        
            for y = 1:length(DATA)
                predicted_matrix(x,y) = Dist(coord_all(x,:), coord_all(y,:)); 
            end;
        end;
    end