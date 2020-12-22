from flask import Flask, jsonify, abort, render_template
import os

templates = f"{os.path.dirname(os.path.abspath(__file__))}/templates"
static = f"{os.path.dirname(os.path.abspath(__file__))}/static"
print(static)
app = Flask(__name__, static_folder = static,  template_folder=templates)

@app.route("/")
def welcome():
    #print(__file__)
    return render_template('submit.html')