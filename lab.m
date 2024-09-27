% Specify the path to the Excel file
excelFilePath = 'C:\Users\samue\OneDrive\Desktop\2023_Fall\Lab II reports Naol Samuel\Lab5 Data\Lab#5_Tensile data.xlsx';

% Read data from Excel file
[num, ~, ~] = xlsread(excelFilePath);

% Extract data from the columns
strain = num(:, 2) * (20 / 2.11) * (1e-6);  % Assuming the strain data is in the second column
stress = (num(:, 3) * 10) ./ (12.42 * 1.93);  % Assuming the stress data is in the third column

% Plotting
figure;
plot(strain, stress, 'DisplayName', 'Engineering Stress Vs Strain');
xlabel('Strain');
ylabel('Stress (MPa)');
grid on;
hold on;

% Plot two percent offset line
x = 0.002:0.001:0.006;
y = 57156.94 * (x - 0.002) ;
plot(x, y, 'k--', 'DisplayName', '0.2% Offset Line');

% Legend
legend('show');