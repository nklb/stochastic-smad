function out = getRootDir()
%%-------------------------------------------------
% Returns root directory
% No parameters required
%%-------------------------------------------------
info = what('burstDetection');
if length(info) > 1
    info = info(1);
end
out = info.path(1:end-15);

end

