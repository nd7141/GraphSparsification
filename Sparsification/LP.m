clear all; 

%    fmae = strcat('Flickr/mae01MP.dat');
%    dlmwrite(fmae, []);

for i=1:1:10
    % undirected case
    fA = strcat('Flickr/LP/A', num2str(i), '.dat');
    fb = strcat('Flickr/LP/d', num2str(i), '.dat');
    fD = strcat('Flickr/LP/D', num2str(i), '.dat');

    % directed case
%     fAin = strcat('Flixster/LP/Ain', num2str(i), '.dat');
%     fAout = strcat('Flixster/LP/Aout', num2str(i), '.dat');
%     fdin = strcat('Flixster/LP/din', num2str(i), '.dat');
%     fdout = strcat('Flixster/LP/dout', num2str(i), '.dat');
%     fD = strcat('Flixster/LP/D', num2str(i), '.txt');
    
%     undirected case
    A = load(fA);
    b = load(fb);
    
    % directed case
%     Ain = load(fAin);
%     Aout = load(fAout);
%     din = load(fdin);
%     dout = load(fdout);

    D = load(fD);
%     undirected case
    n = D(1);
    s = D(2);
    A1 = sparse(A(:,1), A(:,2), A(:,3));
    d = sparse(b(:,1), b(:,2), b(:,3), n, 1);

%     directed case
%     
%     n = D(1); % total # of nodes in non-sparsified graph
%     sp_n = D(2);
%     sp_m = D(3);
%     s = D(4); % sum of discreapncies of sparsified nodes
%     
%     Ai = sparse(Ain(:,1), Ain(:,2), Ain(:,3), sp_n, sp_m);
%     Ao = sparse(Ain(:,1), Ain(:,2), Ain(:,3), sp_n, sp_m);
%     di = sparse(din(:,1), din(:,2), din(:,3), sp_n, 1);
%     do = sparse(dout(:,1), dout(:,2), dout(:,3), sp_n, 1);
    
%     d = spconvert(b); % expected degrees of non-sparsified nodes

    m = size(A1,2);

    f = -ones(1,m);
%     A1 = [Ai; Ao];
%     d = [di; do];

%     options = optimset('LargeScale','off','Simplex','on');
    tic
    [x, fval] = linprog(f, A1, d, [], [], zeros(m,1), ones(m,1), zeros(m,1));
    toc 

    a = sum(abs(d - A1*x)) % sum of discreapncies of non-sparsified nodes
    s
    n
    mae = (s + a)/n

    fprob = strcat('Flickr/LP/x', num2str(i), '.dat');
    dlmwrite(fprob, x);
    
%     dlmwrite(fmae, [i/100, mae], '-append');
end