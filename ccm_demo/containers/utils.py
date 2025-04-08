import os

class ContainerRunner:
    def __init__(self, container):
        self.container = container

    def run(self, command, bind_paths=None, gpu=True, cleanenv=True):
        pass

