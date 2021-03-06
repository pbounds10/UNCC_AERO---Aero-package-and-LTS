function tireCoefficients = importTires2(workbookFile, sheetName, dataLines)
%IMPORTFILE Import data from a spreadsheet
%  TIRECOEFFICIENTS = IMPORTFILE(FILE) reads data from the first
%  worksheet in the Microsoft Excel spreadsheet file named FILE.
%  Returns the data as a table.
%
%  TIRECOEFFICIENTS = IMPORTFILE(FILE, SHEET) reads from the specified
%  worksheet.
%
%  TIRECOEFFICIENTS = IMPORTFILE(FILE, SHEET, DATALINES) reads from the
%  specified worksheet for the specified row interval(s). Specify
%  DATALINES as a positive scalar integer or a N-by-2 array of positive
%  scalar integers for dis-contiguous row intervals.
%
%  Example:
%  tireCos = importTires2("C:\Users\Charles\Desktop\College Spring 2020\SD 2\Lap Time Sim\Lap Time 3-19-2020\Lap Time 3-19-2020\Lap Time Simulation\Tires\tireCoefficients.xlsx", "Sheet1", [2, 4]);
%
%  See also READTABLE.
%
% Auto-generated by MATLAB on 10-Apr-2020 14:07:07

%% Input handling

% If no sheet is specified, read first sheet
if nargin == 1 || isempty(sheetName)
    sheetName = 1;
end

% If row start and end points are not specified, define defaults
if nargin <= 2
    dataLines = [2, 4];
end

%% Setup the Import Options
opts = spreadsheetImportOptions("NumVariables", 34);

% Specify sheet and range
opts.Sheet = sheetName;
opts.DataRange = "A" + dataLines(1, 1) + ":AH" + dataLines(1, 2);

% Specify column names and types
opts.VariableNames = ["VarName1", "a1", "a2", "a3", "a4", "a5", "a6", "a7", "a8", "a9", "a10", "a11", "a12", "a13", "shxa", "bxa", "exa", "rcx1", "rcx2", "rcx3", "rcx4", "rcx5", "rcx6", "rcx7","shyk","byk","eyk","rcy1","rcy2","rcy3","rcy4","rcy5","rcy6","rcy7"];
opts.SelectedVariableNames = ["VarName1", "a1", "a2", "a3", "a4", "a5", "a6", "a7", "a8", "a9", "a10", "a11", "a12", "a13", "shxa", "bxa", "exa", "rcx1", "rcx2", "rcx3", "rcx4", "rcx5", "rcx6", "rcx7","shyk","byk","eyk","rcy1","rcy2","rcy3","rcy4","rcy5","rcy6","rcy7"];
opts.VariableTypes = ["string", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];
opts = setvaropts(opts, 1, "WhitespaceRule", "preserve");
opts = setvaropts(opts, 1, "EmptyFieldRule", "auto");

% Import the data
tireCoefficients = readtable(workbookFile, opts, "UseExcel", false);

for idx = 2:size(dataLines, 1)
    opts.DataRange = "A" + dataLines(idx, 1) + ":X" + dataLines(idx, 2);
    tb = readtable(workbookFile, opts, "UseExcel", false);
    tireCoefficients = [tireCoefficients; tb]; %#ok<AGROW>
end

end