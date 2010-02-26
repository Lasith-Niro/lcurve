# Marty Fuhry
# February 20, 2010
# LaTeX File Generator

import sys

# Collect information from the terminal
if len(sys.argv) == 4:
    course     = sys.argv[1]
    assignment = sys.argv[2]
    problems   = sys.argv[3]
        
# Otherwise, ask for input
else:
    course     = raw_input("Course: ")
    assignment = raw_input("Assignment: ")
    problems   = raw_input("Problems: ")

# Generate the LaTeX File using the appropriate parameters

latex = "\
\\documentclass[11pt]{article}\n\
\\usepackage{listings}\n\
\\usepackage[fleqn]{amsmath}\n\
\\usepackage{graphicx}\n\
\\begin{document}\n\
% Start your text\n\
\\newcommand{\makehomework}[2]%\n\
{\\begin{center}%\n\
	\\Huge #1\\\\%\n\
	\\Large #2\\\\%\n\
	Marty Fuhry\\\\%\n\
	\\today%\n\
\\end{center}}\n\
\n\
\\makehomework{"+course+"}{Homework Assignment "+assignment+"}\n\
\n\
\\lstset{language=Matlab,numbers=left,frame=single,breaklines=true,morecomment=[l]{//}}\n\n\
\n\
% Import Program\n\
%\\lstinputlisting{problem.m}\n\
\n\
% Import Graph\n\
%\\begin{center}\n\
%\\includegraphics[scale=0.5]{problem_17_1_graph.png}\n\
%\\end{center}\n\
\n"

for prob in problems:
    try:
        if type(int(prob)) is type(int()):
            latex += "\\section*{Exercise "+prob+"}\n\n"
    except ValueError:
        prob        

latex += "% Stop your text\n\
\\end{document}\n"

# Write to file

file = open(course+"_"+assignment+".tex",'w')
file.write(latex)

# Usage
