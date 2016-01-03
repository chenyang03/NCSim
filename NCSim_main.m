function NCSim_main()
%
% Simulator for Decentralized Network Coordinate Algorithms (NCSim) 
%
% 
% Version 1.1.0
% Updated on Jan. 3, 2016
% 
% Copyright (C) <2011-2016> by Yang Chen, Fudan University (chenyang@fudan.edu.cn)
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in
% all copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
% THE SOFTWARE.
%

%%%%%%%%%%%%%%%%% Code %%%%%%%%%%%%%%%%% 


% raw distance matrix %

clear
load('data_matrix.mat');
% PL: 169 * 169 PlanetLa data set
% Toread: 355 * 355 PlanetLab data set (collected in Mar.-Apr. 2010, 
%   used in our ACM ReArch'10 paper 'Taming the Triangle Inequality Violations with Network Coordinate System on Real Internet',
%   http://code.google.com/p/toread/) 
% kingmatrix: 1740 * 1740 King data set (http://pdos.csail.mit.edu/p2psim/kingdata/)
%   used in many Network Coordinate papers

% DATA = PL; % PlanetLab data set (small)
DATA = Toread; % PlanetLab data set (big)
% DATA = king_matrix; % King data set


% parameters of NCSim
default_dimension = 8; % config your NC dimension here
max_round = 3; % config your round(s) of simulation here
re_cdf_on = 1; % display the CDF of Relative Error (RE)
vivaldi_option = 0;
% selecting Vivaldi branch, 
%   0 - Vivaldi (basic), original Vivaldi
%   1 - Vivaldi (height),  original Vivaldi
%   2 - Vivaldi (TIV aware), used in "Towards Network Triangle Inequality
%   Violation Aware Distributed Systems" 
%   (Proc. of ACM IMC, 2007).
ides_option = 0;
% selecting IDES branch:
%   0 - IDES (nonnegative), ensuring all predicted distances to be nonnegative, 
%   used in "Phoenix: A Weight-based Network Coordinate System Using Matrix Factorization" 
%   (IEEE Transactions on Network and Service Management, 2011, Vol. 8, Issue 4)
%   1 - IDES (SVD)
%   2 - IDES(NMF)

fprintf('\nNCSim (%d Nodes)\n', length(DATA));


if (re_cdf_on == 1)
    %CDF of relative error (RE)

    fprintf('\nCDF of relative error (RE): \n\n', length(DATA));
    total_re_phoenix = []; total_rank_accuracy_phoenix = [];
    total_re_vivaldi = []; total_rank_accuracy_vivaldi = [];
    total_re_dmf = []; total_rank_accuracy_dmf = [];
    total_re_ides = []; total_rank_accuracy_ides = [];
    
    N = length(DATA);

    for round = 1:max_round
        % NC: Phoenix      
        [out_all, in_all, fpre_newhost, fpre_flashcrowd] = NCS_phoenix(DATA, default_dimension, length(DATA), 32, 5, 0, []);        
        predicted_matrix = out_all * in_all; real_matrix = DATA; rerr = relative_error(predicted_matrix, DATA);        
        output_re = store_re(rerr', 1, 1000); total_re_phoenix = [total_re_phoenix; output_re];        
        fprintf('Phoenix: %.2f ', NPRE(rerr));
        phoenix_rank_accuracy = rank_accuracy(predicted_matrix, DATA);        
        total_rank_accuracy_phoenix = [total_rank_accuracy_phoenix; phoenix_rank_accuracy];       
        
        % NC: Vivaldi
   
        % 0- basic Vivaldi, 1- Vivaldi (height), 2-Vivaldi TIV aware
        predicted_matrix = NCS_vivaldi_all(DATA, default_dimension, length(DATA), 32, vivaldi_option); 
        rerr = relative_error(predicted_matrix, DATA);        
        output_re = store_re(rerr', 1, 1000); total_re_vivaldi = [total_re_vivaldi; output_re];
        fprintf('Vivaldi: %.2f ', NPRE(rerr));           
        [vivaldi_rank_accuracy] = rank_accuracy(predicted_matrix, DATA);
        total_rank_accuracy_vivaldi = [total_rank_accuracy_vivaldi; vivaldi_rank_accuracy];
               

        % NC: DMF
        %[out_all, in_all, fperr] = NCS_DMF(DATA, default_dimension, 32, 100, 50, 0);
        [out_all, in_all] = NCS_DMFSGD(DATA, 32);
        predicted_matrix = out_all * in_all; 
        rerr = relative_error(predicted_matrix, DATA);
        %median(rerr)
        output_re = store_re(rerr', 1, 1000);  total_re_dmf = [total_re_dmf; output_re];
        fprintf('DMFSGD: %.2f ', NPRE(rerr));
        [dmf_rank_accuracy] = rank_accuracy(predicted_matrix, DATA);
        total_rank_accuracy_dmf = [total_rank_accuracy_dmf; dmf_rank_accuracy];
%         
        % IDES
        tmp = randperm(N);         
        landmarks = tmp(1:32); % choose random L nodes as landmarks
        hosts = tmp(32+1:N);  % the rest of them are orinary hosts
        D_landmark = DATA(landmarks, landmarks);
        D_host2landmark = DATA(hosts, landmarks);
        
        [out_l, in_l, out_h, in_h] = NCS_IDES_all(D_landmark, D_host2landmark, default_dimension, ides_option);
        predicted_matrix = out_h*in_h;
        real_matrix = DATA(hosts, hosts);
        rerr=relative_error(predicted_matrix, real_matrix);        
        npre_ides = NPRE(rerr);
        output_re = store_re(rerr', 1, 1000); total_re_ides = [total_re_ides; output_re];
        fprintf('IDES: %.2f ', NPRE(rerr));
        [ides_rank_accuracy] = rank_accuracy(predicted_matrix, real_matrix);
        total_rank_accuracy_ides = [total_rank_accuracy_ides; ides_rank_accuracy];        
        
        fprintf('\n');
    end
    
    figure;
    h1 = plot(0:1/1000:1-1/1000, mean(total_re_phoenix), 'b--');set(h1, 'LineWidth', 2);hold on;
    h2 = plot(0:1/1000:1-1/1000, mean(total_re_vivaldi), 'g:');set(h2, 'LineWidth', 2);hold on;
    h3 = plot(0:1/1000:1-1/1000, mean(total_re_dmf), 'k-');set(h3, 'LineWidth', 2);hold on;
    h4 = plot(0:1/1000:1-1/1000, mean(total_re_ides), 'r-.');set(h4, 'LineWidth', 2);hold on;
    h0 = plot(0:1, [0.9 0.9], 'r:');hold on;
    
    xlabel('Relative Error', 'FontSize', 16);ylabel('Cumulative Distribution Function', 'FontSize', 16);axis([0 1 0 1]);
    h5=legend('Phoenix', 'Vivaldi', 'DMFSGD', 'IDES', 'Location', 'SouthEast');set(h5, 'FontSize', 16);
    RE_filename = 'NCSim_RECDF_';
    tmp_size = length(DATA);
    RE_filename = strcat(RE_filename, num2str(tmp_size));
    saveas(gcf, RE_filename, 'eps');


    
    figure;
    percentage_vec = [0.01:0.01:0.1 0.2:0.1:1];
    h1 = semilogx(percentage_vec, mean(total_rank_accuracy_phoenix), 'b--');set(h1, 'LineWidth', 2);hold on;
    h2 = semilogx(percentage_vec, mean(total_rank_accuracy_vivaldi), 'g:');set(h2, 'LineWidth', 2);hold on;
    h3 = semilogx(percentage_vec, mean(total_rank_accuracy_dmf), 'k-');set(h3, 'LineWidth', 2);hold on;
    h4 = semilogx(percentage_vec, mean(total_rank_accuracy_ides), 'r-.');set(h4, 'LineWidth', 2);hold on;
    xlabel('Fraction of Shortest Paths to Predict (Log Scale)', 'FontSize', 16);ylabel('Cumulative Distribution Function', 'FontSize', 16);
    axis([1/100*0.9 1 0 1]);
    h5=legend('Phoenix', 'Vivaldi', 'DMFSGD', 'IDES',  'Location', 'SouthEast');set(h5, 'FontSize', 16);
    RE_filename = 'Ranking_Accuracy_';
    tmp_size = length(DATA);
    RE_filename = strcat(RE_filename, num2str(tmp_size));
    saveas(gcf, RE_filename, 'eps');

    
    phoenix_fpre = seek_percentage(mean(total_re_phoenix), 0.5);phoenix_npre = seek_percentage(mean(total_re_phoenix), 0.9);
    vivaldi_fpre = seek_percentage(mean(total_re_vivaldi), 0.5);vivaldi_npre = seek_percentage(mean(total_re_vivaldi), 0.9);
    dmf_fpre = seek_percentage(mean(total_re_dmf), 0.5);dmf_npre = seek_percentage(mean(total_re_dmf), 0.9);
    ides_fpre = seek_percentage(mean(total_re_ides), 0.5);ides_npre = seek_percentage(mean(total_re_ides), 0.9);
    fprintf('\nVivaldi|IDES|DMFSGD|Phoenix');
    fprintf('\n50th Percentile RE: %.3f %.3f %.3f %.3f\n', vivaldi_fpre, ides_fpre, dmf_fpre, phoenix_fpre);
    fprintf('90th Percentile RE: %.3f %.3f %.3f %.3f\n', vivaldi_npre, ides_npre, dmf_npre, phoenix_npre);
    
    
end
