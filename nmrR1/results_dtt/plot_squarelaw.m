% define which carbons to plot on a square-law plot
carbons = ["C29","C210","C211","C39","C310","C311"];

% read the carbons' order parameters and relaxation rates
% make sure to specify the correct paths to the files
data = zeros(length(carbons),2);
bexp = [];
for i=1:length(carbons)
    scd = load(strcat('../scd_',carbons(i),'.dat'));
    R1 = load(strcat('R1_',carbons(i),'.dat'));
    data(i,:) = [scd^2 R1];
end

% plot the relaxation rate as a function of the squared order parameter
figure
hold on
plot(data(:,1),data(:,2),'o','LineWidth',1);
xlabel('|S_{CD}|^2');
ylabel('R_{1Z} [1/s]');
