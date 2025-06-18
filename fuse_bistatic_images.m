function I_fused = fuse_bistatic_images(I_CS, I, method)
% Performs coherent or magnitude fusion of radar pair images
%FUSE_BISTATIC_IMAGES Coherently fuse monostatic and bistatic images.
%   I_fused = fuse_bistatic_images(I_CS, METHOD) combines the images stored
%   in cell array I_CS. If METHOD is 'coherent' (default) the images are
%   summed coherently using a global phase alignment between bistatic pairs.
%   Otherwise a magnitude-only summation is performed.
%
%   I_CS{n,m} should contain the range-angle map for transmitter n and
%   receiver m after coarse synchronization.
%
%   The routine first sums the monostatic images (n==m). For each bistatic
%   pair {n,m} with n<m, the corresponding image {m,n} is phase aligned by
%   estimating a single complex phase factor and then the two contributions
%   are added. All results are accumulated into I_fused which is finally
%   normalised.
%
%   Example:
%       I_fused = fuse_bistatic_images(I_CS, 'coherent');
%
%   See also GENERATE_BISTATIC_IMAGES

if nargin < 2
    method = 'coherent';
end

N = size(I_CS,1);
I_fused = zeros(size(I_CS{1,1}));

% --- Sum monostatic images ---
for n = 1:N
    I_fused = I_fused + I{n,n};
end

% --- Fuse bistatic pairs ---
for n = 1:N
    for m = n+1:N
        I_nm = I_CS{n,m};
        I_mn = I_CS{m,n};

        if strcmpi(method,'coherent')
            % Estimate global phase difference between the two images
            phase_diff = angle(sum(I_nm(:) .* conj(I_mn(:))));
            I_mn_aligned = I_mn .* exp(-1j*phase_diff);
            I_pair = I_nm + I_mn_aligned;
        else
            % Non-coherent (magnitude) summation
            I_pair = abs(I_nm) + abs(I_mn);
        end

        I_fused = I_fused + I_pair;
    end
end

% Normalise output
max_val = max(abs(I_fused(:)));
if max_val > 0
    I_fused = I_fused ./ max_val;
end
end
