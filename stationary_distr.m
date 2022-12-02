% This code illustrates how to find the Invariant (Stationary-Ergodic) Distribution of a Markov Chain
% Based on Chris Edmond's notes
% Panagiotis Veneris, U of Liverpool
% 23/1/2021

% The Markov chain (x,P,pi) is given by
% x  : state space (states of the chain, ex 2-state markov chain: low-high income)
% P  : probability transition matrix
% pi : distribution of the variable

% pi_bar = P' * pi_bar ( or pi_bar' = pi_bar'*P as in Sargent-Ljungqvist)
% pi_bar = (stationary) probability distribution


clear all; close all; clc;


% Transition matrix of Markov Chain

P = [0.7 0.2 0.1; 0 0.5 0.5; 0 0.9 0.1];


% Step1:  Compute matrix of eigenvalues and eigenvectors of (the transposed transition matrix) P': remember the transpose!!!! 

[V,D] = eig(P'); % V:matrix of eigenvectors (of transition matrix)
                 % D:diagonal matrix containing the eigenvalues of (transition matrix) P' on the main diagonal
                 
                 
                 
% Step 1.1: Verify that P'V=VD (T=TT)   
T = P'*V;
TT = V*D;


% Step 2: Find which column contains a unit eigenvalue

[x,i] = min(abs(diag(D)-1));  % i=2 -> 2nd column of D contains unitary eigenvalue 

% find the column i of the diagonal matrix D for which
% x=abs(diag(D)-1)=0, that is, for which the eigenvalue is 1


% Step 3: Pick the column from matrix of eigenvectors V (of the transition matrix) associated with
% the unit eigenvalue (2nd column of V)

v = V(:,i);


% Step4: normalize eigenvector found above to sum to 1

X = v/sum(v); % that's the stationary distribution (pi_bar)



% The stationary distribution of a Markov chain is given by the eigenvector
% (of the prob. transition matrix) associated to its unitary eigenvalue 
% (the eigenvector is v but after normalization to sum to 1 it's given by X)
