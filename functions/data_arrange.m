function Ybar = data_arrange(data)

M = size(data{1},1);
noSubjs = length(data);
Ybar = cell(1,noSubjs);
for subj = 1:noSubjs
    Ybar{subj}(:,:,1) = zeros(M,M^2);
    for t = 1:size(data{subj},2)-1
        Ybart = zeros(M,M^2);
        st = 1;
        for m = 1:M
            yt_1 = data{subj}(:,t);
            Ybart(m,st:st+M-1) = yt_1';
            st = st + M;
        end
        Ybar{subj}(:,:,t+1) = Ybart;
    end
end
