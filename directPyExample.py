import os
import d3d
from d3dc import *

try:
    import wx
except:
    d3d.messageBox(u"wxPython sample",
                   u"wxPython is not available.\n\nThis sample will not work without it.")
    raise
# Simple pre-transformed triangle.
triangle = (
    (10.0, 400.0, 0.5, 1.0, 0x80ffff00),
    (190.0, 10.0, 0.5, 1.0, 0x80ff00ff),
    (370.0, 400.0, 0.5, 1.0, 0x8000ff00),
)


class MainFrame(wx.Frame):
    def __init__(self, parent, id, title, pos, size=(400, 600)):
        wx.Frame.__init__(self, parent, id, title, pos, size)

        self.mainpanel = wx.Panel(self, -1)

        self.leftwindow = wx.Panel(self.mainpanel, -1, (10, 10), (380, 450))
        self.leftwindow.SetBackgroundColour(wx.WHITE)
        useleft = wx.Button(self.mainpanel, -1, "Rotate", (100, 500))
        self.Bind(wx.EVT_BUTTON, self.useLeft, useleft)
        self.Bind(wx.EVT_PAINT, self.onPaint)
        self.Bind(wx.EVT_IDLE, self.onIdle)

        # Create a device that preserves its backbuffer after
        # present() by using SWAPCOPY
        d3d.createDevice(u"", u"", 0, 0, False, CREATE.SOFTWARE | CREATE.SWAPCOPY,
                         self.leftwindow.GetHandle())

        self.override = 0
        self.media = None
        self.mesh = d3d.StaticMesh(u"xfiles/static.bigship.x")
        center, radius = self.mesh.getBoundingVolume()
        self.scale = 10.0 / radius
        self.rotation = 0.0
        self.render()  # We only have to render once, use present() to show it.

    def render(self):
        d3d.clear()
        d3d.beginScene()
        d3d.setState(RS.LIGHTING, True)
        d3d.setState(RS.AMBIENT, 0x00202020)
        d3d.setState(RS.SPECULARENABLE, True)
        d3d.setLight(0, LIGHT.POINT, 0xffff0000, 50.0, (0.0, 30.0, 0.0))

        d3d.setState(RS.ZWRITEENABLE, False)
        d3d.setState(RS.ALPHABLENDENABLE, True)
        d3d.setState(RS.FVF, FVF.XYZRHW | FVF.DIFFUSE)
        d3d.drawVertices(TYPE.TRIANGLELIST, triangle)
        d3d.setState(RS.ALPHABLENDENABLE, False)
        d3d.setState(RS.ZWRITEENABLE, True)
        d3d.setView((0.0, 28.0, -28.0,), (0.0, 0.0, 0.0))
        d3d.setTransform((0, 0, 0), (0, self.rotation, 0), (self.scale, self.scale, self.scale))
        self.mesh.render()
        d3d.endScene()

    def useLeft(self, evt):
        # self.override = 0
        self.rotation += 0.1
        self.render()
        # self.leftwindow.Refresh()

    def present(self):
        try:
            d3d.present(self.override)
        except:
            # Device is not functional. For example
            # Fast user switching or power saving modes can cause this.
            try:
                # Try to reset.
                d3d.reset()
            except:
                # Failed, try again later.
                return
            # All ok, continue as usual. The backbuffer has
            # been lost however, so we need to render it again.
            self.render()

    def onIdle(self, evt):
        self.present()

    def onPaint(self, evt):
        dc = wx.PaintDC(self)
        # self.present()


class MainApp(wx.App):
    def OnInit(self):
        frame = MainFrame(None, -1, "DirectPython with wxPython", (10, 10), (400, 600))
        frame.CenterOnScreen()
        frame.Show()
        return True


if __name__ == "__main__":
    app = MainApp(False)
    app.MainLoop()

