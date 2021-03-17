
import xml.etree.ElementTree as ET

tree = ET.parse('parameters.xml')
root = tree.getroot()

print root.find('./mouse').text