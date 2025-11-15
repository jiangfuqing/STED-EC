#!/usr/bin/env python
# coding: utf-8

from PIL import Image
import argparse
import re

ap = argparse.ArgumentParser()
ap.add_argument("-in", "--input_image", required=True, help="input_image")
ap.add_argument("-out", "--output_image", required=True, help="output_image")

args = vars(ap.parse_args())

input_image = args["input_image"]
output_image = args["output_image"]


def resize_image(input_image, output_image):
    img = Image.open(input_image)
    # obtain the image original size
    width, height = img.size
    # calculate new size
    if width > height:
        new_height = 1000
        new_width = int(new_height * width / height)
    else:
        new_width = 1000
        new_height = int(new_width * height / width)
    # resize the image
    resized_img = img.resize((new_width, new_height), Image.ANTIALIAS)
    # save image
    resized_img.save(output_image)

resize_image(input_image, output_image)