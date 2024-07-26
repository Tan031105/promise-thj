import json

def urls_to_json(urls):
    url_list = urls.strip().split('\n')
    json_list = [{"url": url.strip(), "alt": ""} for url in url_list if url.strip()]
    return json_list

if __name__ == "__main__":
    print("请输入多个URL，每个URL换行分隔，输入结束后按Ctrl+D（Linux/Mac）或Ctrl+Z（Windows）并回车:")
    
    # 读取用户输入
    input_urls = ""
    try:
        while True:
            line = input()
            input_urls += line + "\n"
    except EOFError:
        pass
    
    # 转换为JSON格式
    json_output = urls_to_json(input_urls)
    
    # 输出到JSON文件
    with open("output_urls.json", "w", encoding="utf-8") as json_file:
        json.dump(json_output, json_file, ensure_ascii=False, indent=4)
    
    print("URL已成功写入到output_urls.json文件中。")
