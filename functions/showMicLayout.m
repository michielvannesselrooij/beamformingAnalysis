function showMicLayout(x, y)
% CREATES FIGURE WITH MICROPHONE POSITIONS AND NUMBERS

figure;

% Plot data
scatter(x, y, 'b.');

% Add number labels
dx = 0.01;
dy = 0.01;

hold on
for i=1:length(x)
    text(x(i)+dx, y(i)+dy, num2str(i))
end

% Format
axis equal
box on
title('Microphone positions')
xlabel('x [m]');
ylabel('y [m]');

% Set limits
marg = 0.1;
xRange = max(x)-min(x);
yRange = max(y)-min(y);
xlim([ min(x) - marg * xRange, max(x) + marg * xRange ]);
ylim([ min(y) - marg * yRange, max(y) + marg * yRange ]);