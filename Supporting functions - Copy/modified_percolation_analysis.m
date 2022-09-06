function [threshold] = modified_percolation_analysis(A)

%%Input
%A= Adjacency matrix

%%Output
%threshold= Absolute threshold value to be applied to A

threshold=realmin('double');
if number_connected_components(A) > 1
    fprintf(2,'Already more than connected components\n');
    threshold = min(A(:));
    return;
end

a = max(A(:));
b = min(A(:));
TOL = 1E-12;
i = 1;
imax = 10000;

condition_iterations = true;
condition_tolerance = true;
while ( condition_iterations && condition_tolerance )
    condition_iterations = i<imax;
    condition_tolerance = abs(b-a)/2 > TOL;
    %fprintf(2,' %g < T < %g %g\n',b,a,b-a);
    c = (a + b)/2;
    if ( (number_connected_components(threshold_absolute(A,c)~=0) -1 == 0) || (b-a)/2 < TOL )
        threshold = c;
    end
    if sign(number_connected_components(threshold_absolute(A,c)~=0)-1) == sign(number_connected_components(threshold_absolute(A,a)~=0)-1)
        a = c;
    else
        b = c;
    end
    i = i + 1;
end

