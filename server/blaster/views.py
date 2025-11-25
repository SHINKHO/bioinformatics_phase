from django.shortcuts import render
from django.http import HttpRequest, HttpResponse

def index(request: HttpRequest) -> HttpResponse:
    """
    Index view for the blaster app
    """
    return HttpResponse("testPage");
