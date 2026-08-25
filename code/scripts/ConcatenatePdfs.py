#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
#
#
######################################################################
# @author      : bjf79 (bjf79@midway2-login1.rcc.local)
# @file        : ConcatenatePdfs
# @created     : Friday Oct 18, 2024 19:52:37 CDT
#
# @description : 
######################################################################

import sys
import PyPDF2
from PyPDF2 import PdfReader, PdfWriter
from reportlab.pdfgen import canvas
from reportlab.lib.pagesizes import A4
from io import BytesIO

def merge_pdf_pages(output_path, input_paths):
    output_writer = PdfWriter()
    
    for pdf_path in input_paths:
        # Open the input PDF
        with open(pdf_path, "rb") as f:
            reader = PdfReader(f)
            for page_num in range(len(reader.pages)):
                # Get the page object
                page = reader.pages[page_num]

                # Add the page to the output
                output_writer.add_page(page)

    # Save the merged PDF to a file
    with open(output_path, "wb") as output_pdf:
        output_writer.write(output_pdf)

def main(output_path, *input_paths):
    if len(input_paths) < 1:
        print("Please provide at least one input PDF file.")
        sys.exit(1)
    
    merge_pdf_pages(output_path, input_paths)

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python ConcatenatePdfs.py <output.pdf> <InputPlot1.pdf> <InputPlot2.pdf> ...")
        sys.exit(1)
    
    output_file = sys.argv[1]
    input_files = sys.argv[2:]
    
    main(output_file, *input_files)
