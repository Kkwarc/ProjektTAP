function [M, MP] = DMCmatrices(S, N, Nu)
    [ny, nu, D] = size(S);
    M=zeros(N*ny, Nu*nu);
    for i=1:N
        M((i-1)*ny+1:(i-1)*ny+ny, 1:nu) = S(:, :,min(i, D)); 
    end
    for i=2:Nu
        M(:, (i-1)*nu+1:(i-1)*nu+nu) = [zeros((i-1)*ny, nu); M(1:(N-i+1)*ny,1:nu)];
    end
    
    MP=zeros(ny*N, nu*(D-1));
    % for i=1:D-1
    %     for j=1:N
    %         MP((j-1)*ny+1:j*ny,(i-1)*nu+1:i*nu) = S(:, :, min(i+j, D)-S(:, :, i));
    %     end
    % end
    for i = 1:D-1
        for j = 1:N
            index = max(1, min(i+j, D));  % Calculate index here
            MP((j-1)*ny+1:j*ny, (i-1)*nu+1:i*nu) = S(:, :, index)-S(:, :, i);
        end
    end
end